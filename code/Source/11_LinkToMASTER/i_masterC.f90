!
      SUBROUTINE i_masterC(i_opt)
!
!     This routine initializes the coupling-related variables
!
      USE Vol_DATA
      USE Zconst1     , ONLY: cplmaster
      USE Zrv_ncell   , ONLY: p3d_cupid,qvol_mas,ncell_fuel_rod,cupid_cell_hts2d, &
                              master_to_assem1,master_to_assem1_cell,              &
                              master_to_assem1_rod,master_to_assem1_cell,          &
                              asm_ni,asm_nz,chn_nx,chn_ny
!                              
      USE zcore      , ONLY: myrank,np
      USE zmpi       , ONLY: jperm
      USE zzone      , ONLY: ncell_fluid
      USE Ztimecon   , ONLY: time
      USE Zqvol      , ONLY: qvol_liq
      USE Zcoord3    , ONLY: volp
!
      USE MASTER4    , ONLY: TTIME,NZ_TH,NXYF,NPINX,TFC,TFS,TCOO,       &
                             DCOO,BCOO,P3D_TH,VOL_TH,PIN3D_TH,ZCB,NCB, &
                             CBNAM,DTTR,PWTH,power_master,PPM,         &
                             ppm_mas,i_flag,NXY_TH,pin3d_th,           &
                             rv_model_master
!DEC$IF defined (master_flag) 
      USE MASTER4    , ONLY: ppct_master,RSTFN_INP,RSTFN_REV,iok
!DEC$ENDIF 
!
      IMPLICIT none
      INCLUDE "master_c.h"
      INCLUDE "master_dll.h"
!
      INTEGER i_opt
      INTEGER i,j,k
      INTEGER iz,im,irod
!
      REAL(8) p3d_cupid0
!
      INTEGER ii,ij,icupid,jcupid,ni,nz,cx,cy
!
      TTIME=time
!
!.....MASTER4 Initialize
!
!.....Allocate TH input argument
      NXY_TH=(177+64)*4
      NZ_TH =26
      NXYF=177
      NPINX=16
      ALLOCATE(TFC(NXY_TH,NZ_TH),TFS(NXY_TH,NZ_TH))
      ALLOCATE(TCOO(NXY_TH,NZ_TH),DCOO(NXY_TH,NZ_TH),BCOO(NXY_TH,NZ_TH))
      ALLOCATE(P3D_TH(NXY_TH,NZ_TH),VOL_TH(NXY_TH,NZ_TH))
      ALLOCATE(PIN3D_TH(NPINX,NPINX,NXYF,NZ_TH))
!
!.....READ Control Rod information
      NCB=11
      ALLOCATE(ZCB(NCB),CBNAM(NCB))
      CBNAM=(/"R1  ", &
              "S   ", &
              "R2  ", &
              "R3  ", &
              "B   ", &
              "B12 ", &
              "R4  ", &
              "P   ", &
              "R5  ", &
              "A   ", &
              "B11 "   /)
      DTTR=0.0d0
      PWTH=power_master
      PPM=ppm_mas !10.1408942623417
      ZCB=381.0d0
!
!.....for initialization
      i_flag=0
!
!.....CUPIDMASTER (MASTER4)
!
      DO j=1,NZ_TH
         DO i=1,NXY_TH
            P3D_TH(i,j)=0.0d0
            VOL_TH(i,j)=0.0d0
            TFC(i,j)   =0.0d0
            TFS(i,j)   =0.0d0
            TCOO(i,j)  =0.0d0
            DCOO(i,j)  =0.0d0      
         ENDDO
      ENDDO
      PIN3D_TH=0.0d0
!
      IF(myrank.eq.0)THEN
!DEC$IF defined (master_flag) 
         CALL MASTER_DLL(I_FLAG, TTIME, DTTR,                        &
                         TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                         PWTH,                                       &
                         NXY_TH, NZ_TH,                              &
                         PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                         IOK, RSTFN_INP)
         WRITE(*,*) 'i_masterC: Loading 3D Kinetics DLL (master.dll)'
!DEC$ENDIF 
      ENDIF
!
      DO i=1,NXY_TH
         DO j=1,NZ_TH
            DCOO(i,j)=735.18d0 !cell%rhom(2)         ! kg/m3
            TCOO(i,j)=296.21d0 !cell%tl(2)-273.15d0  ! C
            TFC(i,j)=TCOO(i,j)+300.0d0 !-273.15d0     ! C
            TFS(i,j)=TCOO(i,j)+  5.0d0   !-273.15d0     ! C
         ENDDO
      ENDDO
!         
!......Call MASTER DLL
! 
      IF(myrank.eq.0)THEN
         i_flag=i_opt
!DEC$IF defined (master_flag) 
         IF(i_flag.eq.1)THEN
            WRITE(*,*) 'i_masterC: Ready for steady state calculation'
            WRITE(97,*) 'i_masterC: Ready for steady state calculation'
            DO WHILE(1)
               CALL MASTER_DLL(1, TTIME, DTTR,                        &
                               TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                               PWTH,                                       &
                               NXY_TH, NZ_TH,                              &
                               PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                               IOK, RSTFN_INP)
               IF(iok.eq.2)EXIT
            ENDDO                
            WRITE(*,*) 'i_masterC: Finish steady state calculation'
            WRITE(97,*) 'i_masterC: Finish steady state calculation'
!
!........Restart in operation mode  
         ELSEIF(i_flag.eq.3)THEN
            RSTFN_REV="MAS_RST.STD.BIN"
            WRITE(*,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            WRITE(97,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            CALL MASTER_DLL(3, TTIME, DTTR,                        &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_REV)
            WRITE(*,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            WRITE(97,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            i_opt=2
!                            
!........Restart in initialization mode 
         ELSEIF(i_flag.eq.4)THEN
            RSTFN_REV="MAS_RST.STD.BIN"
            WRITE(*,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            WRITE(97,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            CALL MASTER_DLL(3, TTIME, DTTR,                        &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_REV)    
            WRITE(*,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            WRITE(97,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            i_opt=1
         ELSE
            WRITE(*,"(a,1i3)")'i_opt should be 1,3,4 in i_masterC !!!',i_opt
            WRITE(97,"(a,1i3)")'i_opt should be 1,3,4 in i_masterC !!!',i_opt              
            PAUSE
            STOP                     
         ENDIF    
!DEC$ENDIF 
      ENDIF !myrank.eq.0      
      CALL barrier_mpi
!
!.....Assembly-scale
!     p3d_cupid: W/cm3 (volumetric)--> Need to convert to total power (W)   
      IF(cplmaster.eq.1)THEN
!
!........Communicate P3D_TH, vol_TH
!
         IF(np.gt.1)THEN
            DO i=1,NZ_TH
               CALL allreducei_r(P3D_TH(1,i),NXY_TH)
               CALL allreducei_r(vol_TH(1,i),NXY_TH)
            ENDDO
         ENDIF

         p3d_cupid(:)=0.0d0
         qvol_mas(:)=0.0d0
         master_to_assem1_rod=0
         master_to_assem1_cell=0
         DO iz=1,NZ_TH
            DO im=1,NXY_TH !2x2 subchannel
               icupid=master_to_assem1(im,iz)
               do irod=1,ncell_fuel_rod ! 1 assembly
                  k=cupid_cell_hts2d(irod)
                  if(jperm(k).eq.icupid) then
                     master_to_assem1_rod(im,iz)=irod
                     master_to_assem1_cell(im,iz)=k
                     p3d_cupid(irod)=p3d_cupid(irod)+P3D_TH(im,iz)*1.e6*vol_th(im,iz)/1.e6
                     qvol_mas(k)=qvol_mas(k)+P3D_TH(im,iz)*vol_th(im,iz)
                  endif
               enddo
            ENDDO
         ENDDO
!         
!........RV OFF
         IF(rv_model_master.eq.0)THEN
            qvol_liq=0.0d0
            p3d_cupid0=0.d0
            do i=1,ncell_fluid
                qvol_liq(i)=qvol_mas(i)/volp(i)
                p3d_cupid0=p3d_cupid0+qvol_liq(i)*volp(i)
            enddo
            IF(np.gt.1) CALL allreducei_r1(p3d_cupid0)
            WRITE(*,*) 'Total power in i_masterC is', p3d_cupid0
         ENDIF   
!      
!.....rod-scale
!     PIN3D_TH: W (Total power per each rod)
      ELSEIF(cplmaster.eq.2)THEN
         IF(np.gt.1)THEN
            DO i=1,NZ_TH
               CALL allreducei_r(PIN3D_TH(1,1,1,i),NPINX*NPINX*NXYF)
            ENDDO
         ENDIF

         p3d_cupid0=0.0d0
         DO iz=1,NZ_TH
            DO im=1,NXYF
               DO ij=1,NPINX
                  DO ii=1,NPINX
                     p3d_cupid0=p3d_cupid0+PIN3D_TH(ii,ij,im,iz)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF(myrank.eq.0)write(97,*)'total power from MASTER, PIN3D_TH is',p3d_cupid0
         IF(myrank.eq.0)write(* ,*)'total power from MASTER, PIN3D_TH is',p3d_cupid0

         ! PIN3D_TH: W (total power)
         p3d_cupid0=0.0d0

         ! RV ON
         IF(rv_model_master.eq.1)THEN
            DO irod=1,ncell_fuel_rod
               icupid=cupid_cell_hts2d(irod)
               jcupid=jperm(icupid)
               ni=asm_ni(jcupid)
               nz=asm_nz(jcupid)
               cx=chn_nx(jcupid)
               cy=chn_ny(jcupid)
               p3d_cupid(irod)=PIN3D_TH(cx,cy,ni,nz)
               p3d_cupid0=p3d_cupid0+p3d_cupid(irod)
            ENDDO   
         ! RV OFF
         ELSE
            !CALL p3d_to_qvol
            !DO i=1,ncell_fluid
            !   p3d_cupid_all=p3d_cupid_all+qvol_liq(i)*volp(i)
            !ENDDO
         ENDIF     
            
         IF(np.gt.1) CALL allreducei_r1(p3d_cupid0)
         IF(myrank.eq.0)THEN
            WRITE(* ,*) 'total power from Fuel rod, p3d_cupid is',p3d_cupid0
            WRITE(97,*) 'total power from Fuel rod, p3d_cupid is',p3d_cupid0
            WRITE(* ,*) 'Succeed in MASTER-to-Rod in i_masterC'
            WRITE(97,*) 'Succeed in MASTER-to-Rod in i_masterC'
         ENDIF
!         
      ENDIF !cplmaster==1,2
      CALL barrier_mpi
          
      IF(myrank.eq.0)WRITE(* ,*) 'i_masterC: Finish i_masterC'
      IF(myrank.eq.0)WRITE(97,*) 'i_masterC: Finish i_masterC'
!
      END SUBROUTINE i_masterC
